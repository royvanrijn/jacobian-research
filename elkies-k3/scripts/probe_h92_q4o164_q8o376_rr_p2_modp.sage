#!/usr/bin/env sage -python
"""Search the q8/orbit376 RR pencil in the post-collision projective plane.

The exact q8 horizontal H has P.O=4 and compact coordinates

    x(H)=X/Z^2, y(H)=Y/Z^3, deg(Z)=4.

For the selected degree-two neighbour the fibre twist is zero. A global chord
function may therefore be written

    f = AA/Z^2 + (BB/Z) * m,
    m = (y-y(H))/(x-x(H)),
    deg(AA)<=8, deg(BB)<=2.

Smooth collision regularity is the congruence

    AA*X == BB*Y (mod Z^2).

It has rank eight. Modulo the constant function, every post-collision
candidate is consequently determined by a projective quadratic
BB=b0+b1*t+b2*t^2. This script enumerates P^2(GF(p)), forms the exact chord
square class, and retains directions for which the residual vertical square
content leaves a quartic elliptic K3 Jacobian. A semistable 4A1 fingerprint
is recorded as a strong selector but is not assumed a priori.

This is a finite-field discovery probe. It does not certify the QQ pencil.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, gcd, inverse_mod, is_prime


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
HORIZONTAL = LOCAL / "q4o164-q8o376-horizontal-crt-qq.json"
COST = GENERATED / "elkies-k3-h3-q4o164-c8-q8d2-cap10000-growth-equation-cost.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=41)
parser.add_argument("--output", type=Path)
parser.add_argument(
    "--stop-after",
    type=int,
    help="optional processed-candidate cap for smoke tests",
)
args = parser.parse_args()

p = ZZ(args.prime)
if not is_prime(p) or p in (2, 3):
    raise SystemExit("--prime must be a prime other than 2 or 3")
if args.stop_after is not None and args.stop_after <= 0:
    raise SystemExit("--stop-after must be positive")
OUTPUT = (
    args.output
    if args.output and args.output.is_absolute()
    else ROOT / args.output
    if args.output
    else LOCAL / f"q4o164-q8o376-rr-p2-scan-mod{p}.json"
)
INPUTS = (MODEL, HORIZONTAL, COST)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def log(stage, **fields):
    suffix = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"Q8O376P2|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{suffix}" if suffix else ""),
        flush=True,
    )


def reduce_q(value, field):
    value = QQ(value)
    if value.denominator() % p == 0:
        raise ArithmeticError(f"bad denominator modulo {p}")
    return field(value.numerator()) / field(value.denominator())


def polynomial_from_record(record, key, ring, field):
    return ring([reduce_q(value, field) for value in record[key]])


def exact_power_root(polynomial, exponent):
    """Return the monic exact exponent-th root in a finite polynomial ring."""
    polynomial = polynomial.monic()
    answer = polynomial.parent().one()
    for factor, multiplicity in polynomial.factor():
        multiplicity = int(multiplicity)
        if multiplicity % exponent:
            raise ArithmeticError(
                f"polynomial is not an exact {exponent}-th power: multiplicity {multiplicity}"
            )
        answer *= factor.monic() ** (multiplicity // exponent)
    if answer**exponent != polynomial:
        raise ArithmeticError("exact polynomial power-root verification failed")
    return answer.monic()


def projective_quadratics(field, variable):
    """Canonical representatives of P^2(field), first nonzero coefficient one."""
    for b1 in field:
        for b2 in field:
            yield field.one() + b1 * variable + b2 * variable**2, (1, int(b1), int(b2))
    for b2 in field:
        yield variable + b2 * variable**2, (0, 1, int(b2))
    yield variable**2, (0, 0, 1)


def square_root_if_square(polynomial):
    if not polynomial:
        return None
    polynomial = polynomial.monic()
    answer = polynomial.parent().one()
    for factor, multiplicity in polynomial.factor():
        multiplicity = int(multiplicity)
        if multiplicity % 2:
            return None
        answer *= factor.monic() ** (multiplicity // 2)
    return answer.monic() if answer**2 == polynomial else None


def projective_chart(projective):
    return next(index for index, value in enumerate(projective) if value)


model = json.loads(MODEL.read_text())
horizontal = json.loads(HORIZONTAL.read_text())
cost = json.loads(COST.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert horizontal["status"] == "PASS_EXACT_QQ_Q4O164_Q8O376_HORIZONTAL_CRT"
assert cost["status"] == "PASS_EXACT_MARKED_FRONTIER_EQUATION_COST_SCORING"
selected_cost = next(
    record
    for record in cost["retained_candidates"]
    if record["candidate_id"] == {"q": 8, "old_fibre_degree": 2, "orbit_index": 376}
)
assert selected_cost["horizontal"]["P_dot_O"] == 4
assert selected_cost["horizontal"]["fibre_twist"] == 0
assert selected_cost["horizontal"]["vertical_layers"] == 2
assert selected_cost["expected_RR_ambient"] == 12

F = GF(p)
FT = PolynomialRing(F, "t")
t = FT.gen()
FU = PolynomialRing(F, "u")
u = FU.gen()
RT = PolynomialRing(FU, "t")

A = FT([reduce_q(value, F) for value in model["compact_model"]["A_coefficients_low_to_high"]])
B = FT([reduce_q(value, F) for value in model["compact_model"]["B_coefficients_low_to_high"]])
section = horizontal["section"]
x_num = polynomial_from_record(section["x"], "numerator_coefficients_low_to_high", FT, F)
x_den = polynomial_from_record(section["x"], "denominator_coefficients_low_to_high", FT, F)
y_num = polynomial_from_record(section["y"], "numerator_coefficients_low_to_high", FT, F)
y_den = polynomial_from_record(section["y"], "denominator_coefficients_low_to_high", FT, F)
if not x_den or not y_den:
    raise ArithmeticError("horizontal denominator vanished modulo p")
x_den = x_den.monic()
y_den = y_den.monic()
Zx = exact_power_root(x_den, 2)
Zy = exact_power_root(y_den, 3)
if Zx != Zy:
    raise ArithmeticError("x/y denominators do not share one projective Z")
Z = Zx
X = x_num
Y = y_num
if (X.degree(), Z.degree(), Y.degree()) != (12, 4, 18):
    raise ArithmeticError("q8 horizontal degree profile changed modulo p")
if Y**2 != X**3 + A * X * Z**4 + B * Z**6:
    raise ArithmeticError("q8 horizontal misses compact parent modulo p")
if gcd(X, Z) != 1:
    raise ArithmeticError("horizontal pole divisor is not smooth/coprime modulo p")

Z2 = Z**2
X_inverse = inverse_mod(X, Z2)
residue = FT((Y * X_inverse) % Z2)
if (X * residue - Y) % Z2:
    raise ArithmeticError("smooth collision residue failed")


def lift(poly):
    return RT([FU(value) for value in FT(poly).list()])


A_nested, X_nested, Y_nested, Z_nested = map(lift, (A, X, Y, Z))
Z2_nested = Z_nested**2
Z6_nested = Z_nested**6
log(
    "LOAD",
    prime=p,
    parent_degrees=f"{A.degree()},{B.degree()}",
    horizontal_degrees="12,8,18,12",
    post_collision_dimension=4,
)


def coefficients_in_u(poly):
    """Transpose RT=FU[t] into coefficient polynomials in FT for u."""
    if not poly:
        return [FT.zero()]
    u_degree = max(
        (coefficient.degree() for coefficient in poly.list() if coefficient),
        default=0,
    )
    return [
        FT([F(coefficient[index]) for coefficient in poly.list()])
        for index in range(u_degree + 1)
    ]


def content_in_t(poly):
    coefficients = [item for item in coefficients_in_u(poly) if item]
    if not coefficients:
        return FT.zero()
    answer = coefficients[0]
    for item in coefficients[1:]:
        answer = answer.gcd(item)
        if answer.degree() == 0:
            break
    return answer.monic()


def divide_by_t_polynomial(poly, divisor):
    quotient, remainder = poly.quo_rem(lift(divisor))
    if remainder:
        raise ArithmeticError("nested content division left a remainder")
    return quotient


def child_profile(quartic):
    """Return an exact modular fibre profile, including the fibre at infinity."""
    values = list(quartic.list()) + [FU.zero()] * 5
    e, d, c, b, a = values[:5]
    I = 12 * a * e - 3 * b * d + c**2
    J = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
    A_child = FU(-27 * I)
    B_child = FU(-27 * J)
    Delta = FU(-16 * (4 * A_child**3 + 27 * B_child**2))
    if not Delta or A_child.degree() > 8 or B_child.degree() > 12 or Delta.degree() > 24:
        return None

    degree_A = int(A_child.degree())
    degree_B = int(B_child.degree())
    degree_delta = int(Delta.degree())
    infinity_order = 24 - degree_delta
    factors = list(Delta.factor())
    multiplicity_profile = sorted(
        [int(factor.degree()), int(multiplicity)]
        for factor, multiplicity in factors
    )
    finite_semistable = all(
        int(multiplicity) in (1, 2)
        and A_child.gcd(factor) == 1
        and B_child.gcd(factor) == 1
        for factor, multiplicity in factors
    )
    infinity_multiplicative = bool(
        infinity_order == 0
        or (
            infinity_order in (1, 2)
            and degree_A == 8
            and degree_B == 12
        )
    )
    infinity_kodaira = (
        "smooth"
        if infinity_order == 0
        else f"I{infinity_order}"
        if infinity_multiplicative
        else "non-semistable"
    )
    double_degree = sum(
        int(factor.degree())
        for factor, multiplicity in factors
        if int(multiplicity) == 2
    )
    simple_degree = sum(
        int(factor.degree())
        for factor, multiplicity in factors
        if int(multiplicity) == 1
    )
    geometric_root_rank = double_degree + (1 if infinity_order == 2 else 0)
    geometric_I1_count = simple_degree + (1 if infinity_order == 1 else 0)
    semistable = finite_semistable and infinity_multiplicative
    semistable_4A1 = bool(
        semistable
        and infinity_order in (0, 1, 2)
        and geometric_root_rank == 4
        and geometric_I1_count == 16
    )
    return {
        "A_coefficients_low_to_high": [int(value) for value in A_child.list()],
        "B_coefficients_low_to_high": [int(value) for value in B_child.list()],
        "degrees_A_B_Delta": [degree_A, degree_B, degree_delta],
        "finite_discriminant_factor_degrees_and_multiplicities": multiplicity_profile,
        "infinity_delta_order": infinity_order,
        "infinity_kodaira": infinity_kodaira,
        "semistable": semistable,
        "semistable_4A1_fingerprint": semistable_4A1,
        "finite_double_support_degree": double_degree,
        "finite_simple_support_degree": simple_degree,
        "geometric_root_rank": geometric_root_rank,
        "geometric_I1_count": geometric_I1_count,
        "euler_number": degree_delta + infinity_order,
    }


candidates = []
processed = 0
collision_failures = 0
z6_failures = 0
quartic_failures = 0
unclassified_child_profiles = 0
for BB, projective in projective_quadratics(F, t):
    if args.stop_after is not None and processed >= args.stop_after:
        break
    processed += 1
    AA = FT((BB * residue) % Z2)
    if (AA * X - BB * Y) % Z2:
        collision_failures += 1
        continue

    AA_nested, BB_nested = map(lift, (AA, BB))
    N = AA_nested - RT(u) * Z2_nested
    Db = -BB_nested
    raw = (
        N**4
        - 6 * X_nested * N**2 * Db**2
        - 8 * Y_nested * N * Db**3
        - 3 * X_nested**2 * Db**4
        - 4 * A_nested * Z_nested**4 * Db**4
    )
    after_collision, remainder = raw.quo_rem(Z6_nested)
    if remainder or after_collision.degree() > 8:
        z6_failures += 1
        continue

    content = content_in_t(after_collision)
    square_root = square_root_if_square(content)
    if square_root is None or content.degree() > 4:
        quartic_failures += 1
        continue
    infinity_square_degree = 4 - int(content.degree())
    if infinity_square_degree < 0 or infinity_square_degree % 2:
        quartic_failures += 1
        continue
    residual_quartic = divide_by_t_polynomial(after_collision, content)
    if residual_quartic.degree() != 4:
        quartic_failures += 1
        continue
    profile = child_profile(residual_quartic)
    if profile is None:
        unclassified_child_profiles += 1

    candidate = {
        "projective_BB_coefficients_low_to_high": list(projective),
        "projective_chart_index": projective_chart(projective),
        "AA_coefficients_low_to_high": [int(value) for value in AA.list()],
        "finite_vertical_square_root_coefficients_low_to_high": [
            int(value) for value in square_root.list()
        ],
        "finite_vertical_square_degree": int(content.degree()),
        "infinity_vertical_square_degree": infinity_square_degree,
        "quartic_coefficients_in_old_t_low_to_high": [
            [int(value[index]) for index in range(value.degree() + 1)]
            if value
            else []
            for value in residual_quartic.list()
        ],
        "child": profile,
    }
    candidates.append(candidate)
    log(
        "CANDIDATE",
        index=len(candidates),
        BB=projective,
        finite_square_degree=content.degree(),
        child_4A1=bool(profile and profile["semistable_4A1_fingerprint"]),
    )

expected_size = int(p**2 + p + 1)
complete = processed == expected_size and args.stop_after is None
strong_candidates = [
    candidate
    for candidate in candidates
    if candidate["child"] and candidate["child"]["semistable_4A1_fingerprint"]
]
selected_candidate = strong_candidates[0] if complete and len(strong_candidates) == 1 else None
status = (
    "CANDIDATE_UNIQUE_4A1_MODP_Q4O164_Q8O376_RR_P2_SCAN"
    if selected_candidate
    else "CANDIDATE_MODP_Q4O164_Q8O376_RR_P2_SCAN"
    if candidates
    else "NO_CANDIDATE_MODP_Q4O164_Q8O376_RR_P2_SCAN"
)

payload = {
    "schema": "elkies-k3.q4o164-q8o376-rr-p2-scan-modp.v2",
    "status": status,
    "prime": int(p),
    "search_space": {
        "kind": "P2 projective BB coefficients after exact collision recurrence",
        "expected_size": expected_size,
        "processed": processed,
        "complete": complete,
        "stop_after": args.stop_after,
    },
    "exact_reduction": {
        "ambient_dimension": 12,
        "smooth_collision_rank": 8,
        "post_collision_dimension": 4,
        "constant_basis_pair": {"AA": "Z^2", "BB": "0"},
        "nonconstant_parameter": "BB=b0+b1*t+b2*t^2 up to projective scale",
        "AA_recurrence": "AA=BB*Y/X mod Z^2",
        "vertical_layers_from_lattice": 2,
        "expected_final_kernel_dimension": 2,
    },
    "diagnostics": {
        "collision_failures": collision_failures,
        "Z6_divisibility_failures": z6_failures,
        "square_content_or_quartic_degree_failures": quartic_failures,
        "unclassified_child_profiles": unclassified_child_profiles,
        "semistable_4A1_candidate_count": len(strong_candidates),
        "candidate_count": len(candidates),
    },
    "selected_candidate": selected_candidate,
    "strong_candidates": strong_candidates,
    "candidates": candidates,
    "proof_boundary": (
        "Finite-field discovery only. The scan uses the exact QQ parent and horizontal reduced "
        "modulo p, solves the smooth collision module exactly, and filters by square vertical "
        "content and a degree-four elliptic-K3 residual. A semistable 4A1 fingerprint includes "
        "smooth, I1, and I2 infinity possibilities and is only a strong selector. A surviving "
        "direction must still be matched across primes, reconstructed over QQ, pointed at "
        "P1229, and certified by exact substitution and marked transport."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
    "runtime_seconds": time.monotonic() - started,
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
log(
    "DONE",
    processed=processed,
    complete=complete,
    candidates=len(candidates),
    strong=len(strong_candidates),
    status=payload["status"],
    output=OUTPUT,
)
