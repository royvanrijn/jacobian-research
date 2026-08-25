#!/usr/bin/env sage -python
"""Resolve the physical q10 RR pencil torically over GF(103).

status: ACTIVE_PROOF
claim: exact resolved RR plane at the pinned prime for the physical q10 edge
outputs: artifacts/local/elkies-k3/q24-2a5-direct-q10-toric-rr-mod103.json

The 15-dimensional coefficient space is cut by the smooth P.O collision to
five dimensions.  In the literal physical divisor one has

  D-O-H = -C0-C3-C4-C5-A0 -2C1-2C2-2C6-C7-A1,

so the two I6 blocks impose thirteen displayed low-order vanishing rows with
combined rank three.  We use the pinned tangent orientations and the chord
through -H, compile the resulting binary quartic, and classify its exact
specialized Jacobian.  The characteristic-103 root rank jumps from the
generic lattice target 9 to 11; that specialization is recorded explicitly
rather than mistaken for the characteristic-zero fibre profile.  No
Groebner basis or elimination is used.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
POOL = LOCAL / "q24-2a5-zero-pole-sections-p103.json"
SHELL = LOCAL / "q24-2a5-zero-pole-shell-match-p103.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q24-2a5-direct-q10-toric-rr-mod103.json",
)
args = parser.parse_args()
output = Path(args.output)
if not output.is_absolute():
    output = ROOT / output

pool = json.loads(POOL.read_text())
shell = json.loads(SHELL.read_text())
assert pool["status"] == "PASS_BOUNDED_MOD103_ZERO_POLE_SECTION_ENUMERATION"
assert shell["status"] == "PASS_EXHAUSTIVE_MOD103_ZERO_POLE_SHELL_EMBEDDINGS_CANONICAL_MARKING"

p = ZZ(pool["prime"])
assert p == 103
F = GF(p)
RT = PolynomialRing(F, "T")
T = RT.gen()
A = RT(pool["surface_mod_103"]["A_coefficients_low_to_high"])
B = RT(pool["surface_mod_103"]["B_coefficients_low_to_high"])
E = EllipticCurve(RT.fraction_field(), [0, 0, 0, A, B])


def point(index):
    record = pool["sections"][index]
    return E(
        RT(record["X_coefficients_low_to_high"]),
        RT(record["Y_coefficients_low_to_high"]),
    )


# Compatible exact-shell mapping 35 binds all five terms in the short word.
mapping = shell["all_complete_mappings_exact_index_to_modular_index"][35]
assert [mapping[str(index)] for index in (45, 1, 24, 14)] == [49, 30, 115, 112]
assert [mapping[str(index)] for index in (0, 2, 3)] == [114, 62, 36]
P230 = point(114) + point(62) + point(36)
H = 2 * point(49) - 2 * point(30) + 3 * point(115) + 2 * P230 + 2 * point(112)
hx, hy = H.xy()
Z = RT(hx.denominator()).monic().sqrt().monic()
X = RT(hx * Z**2)
Y = RT(hy * Z**3)
assert RT(hy.denominator()).monic() == Z**3
assert Y**2 == X**3 + A * X * Z**4 + B * Z**6

# The leading X coefficient cancels modulo 103; this is harmless for the
# modular RR rank and is independently checked by the exact QQ compiler.
assert (X.degree(), Y.degree(), Z.degree()) == (13, 21, 5)


# Smooth collision basis: (Z^2,0), then BB=1,T,T^2,T^3 with
# AA=BB*Y/X modulo Z^2.  X is a unit modulo Z^2.
collision_modulus = Z**2
X_inverse = X.inverse_mod(collision_modulus)
pairs = [(Z**2, RT.zero())]
for degree in range(4):
    BB = T**degree
    AA = RT((BB * Y * X_inverse) % collision_modulus)
    pairs.append((AA, BB))
assert len(pairs) == 5


Kr = FunctionField(F, "r")
r = Kr.gen()
LS = LaurentSeriesRing(Kr, "s", default_prec=11)
s = LS.gen()


def shifted_polynomial(poly, root):
    answer = LS.zero()
    base = LS(Kr(root)) + s
    power = LS.one()
    for coefficient in RT(poly).list():
        answer += Kr(coefficient) * power
        power *= base
    return answer


def shifted_rational(value, root):
    value = RT.fraction_field()(value)
    return shifted_polynomial(value.numerator(), root) / shifted_polynomial(
        value.denominator(), root
    )


def newton_sqrt(value, initial):
    answer = LS(Kr(initial))
    for unused in range(4):
        answer = (answer + value / answer) / 2
    return answer


def functional_rows(values):
    """Expand one identity in F(r) into its exact coefficient rows over F."""
    values = [Kr(value) for value in values]
    nonzero = [value for value in values if value]
    if not nonzero:
        return []
    common = nonzero[0].denominator().parent().one()
    for value in nonzero:
        common = common.lcm(value.denominator())
    numerators = [(value * common).numerator() for value in values]
    maximum_degree = max((value.degree() for value in numerators if value), default=-1)
    answer = []
    for degree in range(maximum_degree + 1):
        row = [
            F(value[degree]) if value and degree <= value.degree() else F.zero()
            for value in numerators
        ]
        if any(row):
            answer.append(row)
    return answer


roots = {
    F(89): {
        "node": F(65),
        "rho": F(35),
        "coefficients": (1, 1, 1, 1, 1),
    },
    F(68): {
        "node": F(90),
        "rho": F(95),
        "coefficients": (1, 2, 2, 2, 1),
    },
}


def toric_rows(root, rho_start, chord_sign, function_sign):
    data = roots[root]
    Aloc = shifted_polynomial(A, root)
    Bloc = shifted_polynomial(B, root)
    center = newton_sqrt(-Aloc / 3, data["node"])
    smoothing = center**3 + Aloc * center + Bloc
    assert int(smoothing.valuation()) == 6
    unit = smoothing / s**6
    Hxloc = shifted_rational(hx, root)
    Hyloc = shifted_rational(hy, root)

    rows = []
    diagnostics = []
    for component, required_order in enumerate(data["coefficients"], start=1):
        aa = LS(r) * s**component
        bb = unit * s**(6 - component) / LS(r)
        yy = (aa + bb) / 2
        ww = (bb - aa) / 2

        rho_series = LS(Kr(rho_start))
        for unused in range(4):
            rho_series -= (
                rho_series**3 - 3 * center * rho_series - ww
            ) / (3 * rho_series**2 - 3 * center)
        xx = center + ww / rho_series
        m = (yy + chord_sign * Hyloc) / (xx - Hxloc)

        values = []
        for AA, BB in pairs:
            values.append(
                shifted_rational(AA / Z**2, root)
                + function_sign * shifted_rational(BB / Z, root) * m
            )
        minimum_valuation = min(int(value.valuation()) for value in values if value)
        diagnostics.append({
            "component": int(component),
            "required_vanishing_order": int(required_order),
            "minimum_basis_valuation": int(minimum_valuation),
        })
        for order in range(required_order):
            rows.extend(functional_rows([value[order] for value in values]))
    return rows, diagnostics


def canonical_kernel(rows):
    matrix_rows = matrix(F, rows)
    return matrix_rows, matrix_rows.right_kernel().basis_matrix().echelon_form()


KU = FunctionField(F, "U")
U = KU.gen()
ST = PolynomialRing(KU, "T")


def kodaira_type(va, vb, vd):
    if vd == 0:
        return None, 0
    if va == 0 or vb == 0:
        return f"I{vd}", max(0, vd - 1)
    if vd == 2:
        return "II", 0
    if vd == 3:
        return "III", 1
    if vd == 4:
        return "IV", 2
    if vd == 6 and va >= 2 and vb >= 3:
        return "I0*", 4
    if vd >= 7 and va == 2 and vb == 3:
        return f"I{vd - 6}*", vd - 2
    if vd == 8 and va >= 3 and vb == 4:
        return "IV*", 6
    if vd == 9 and va == 3 and vb >= 5:
        return "III*", 7
    if vd == 10 and va >= 4 and vb == 5:
        return "II*", 8
    return f"unclassified({va},{vb},{vd})", -1


def minimal_orders(va, vb, vd):
    scale = min(floor(QQ(va) / 4), floor(QQ(vb) / 6))
    return int(va - 4 * scale), int(vb - 6 * scale), int(vd - 12 * scale), int(scale)


def valuation_at(value, prime):
    value = KU(value)
    return int(value.numerator().valuation(prime) - value.denominator().valuation(prime))


def infinity_order(value, weight):
    value = KU(value)
    return int(weight + value.denominator().degree() - value.numerator().degree())


def jacobian_profile(A4, A6):
    discriminant = -16 * (4 * A4**3 + 27 * A6**2)
    prime_set = set()
    for value in (A4, A6, discriminant):
        for polynomial in (KU(value).numerator(), KU(value).denominator()):
            for prime, unused in polynomial.factor():
                prime_set.add(prime.monic())

    places = []
    root_rank = 0
    euler = 0
    ade = []
    for prime in sorted(prime_set, key=lambda q: (q.degree(), str(q))):
        raw = (
            valuation_at(A4, prime),
            valuation_at(A6, prime),
            valuation_at(discriminant, prime),
        )
        va, vb, vd, scale = minimal_orders(*raw)
        if vd == 0:
            continue
        kind, rank = kodaira_type(va, vb, vd)
        if rank > 0:
            root_rank += int(prime.degree()) * rank
        if rank >= 0:
            euler += int(prime.degree()) * vd
        if kind and rank > 0:
            ade.extend([kind] * int(prime.degree()))
        places.append({
            "prime": str(prime),
            "degree": int(prime.degree()),
            "raw_orders_A4_A6_Delta": list(map(int, raw)),
            "minimal_scale": int(scale),
            "minimal_orders_A4_A6_Delta": [va, vb, vd],
            "kodaira": kind,
            "root_rank_per_geometric_place": int(rank),
        })

    raw_infinity = (
        infinity_order(A4, 8),
        infinity_order(A6, 12),
        infinity_order(discriminant, 24),
    )
    va, vb, vd, scale = minimal_orders(*raw_infinity)
    infinity_kind, infinity_rank = kodaira_type(va, vb, vd)
    if infinity_rank > 0:
        root_rank += infinity_rank
        ade.append(infinity_kind)
    if infinity_rank >= 0:
        euler += vd
    return {
        "A4_degree": int(KU(A4).numerator().degree()),
        "A4_denominator_degree": int(KU(A4).denominator().degree()),
        "A6_degree": int(KU(A6).numerator().degree()),
        "A6_denominator_degree": int(KU(A6).denominator().degree()),
        "Delta_degree": int(KU(discriminant).numerator().degree()),
        "Delta_denominator_degree": int(KU(discriminant).denominator().degree()),
        "finite_places": places,
        "infinity": {
            "raw_orders_A4_A6_Delta": list(map(int, raw_infinity)),
            "minimal_scale": int(scale),
            "minimal_orders_A4_A6_Delta": [va, vb, vd],
            "kodaira": infinity_kind,
            "root_rank": int(infinity_rank),
        },
        "euler_number": int(euler),
        "root_rank": int(root_rank),
        "ADE": sorted(ade),
    }


def compile_kernel(kernel):
    pair_rows = []
    for row in kernel.rows():
        AA = sum((row[index] * pairs[index][0] for index in range(len(pairs))), RT.zero())
        BB = sum((row[index] * pairs[index][1] for index in range(len(pairs))), RT.zero())
        pair_rows.append((AA, BB))
    (AA0, BB0), (AA1, BB1) = pair_rows
    N = ST(AA1) - U * ST(AA0)
    D = U * ST(BB0) - ST(BB1)
    raw = (
        N**4 - 6 * ST(X) * N**2 * D**2 - 8 * ST(Y) * N * D**3
        - 3 * ST(X)**2 * D**4 - 4 * ST(A) * ST(Z)**4 * D**4
    )
    quotient, remainder = raw.quo_rem(ST(Z)**6)
    assert not remainder
    factors = quotient.factor()
    odd = ST(factors.unit())
    factor_data = []
    for factor, exponent in factors:
        factor_data.append((int(factor.degree()), int(exponent)))
        if exponent % 2:
            odd *= factor
    record = {
        "kernel_basis": [[int(value) for value in row] for row in kernel.rows()],
        "raw_degree": int(quotient.degree()),
        "raw_factor_degrees_and_exponents": factor_data,
        "odd_degree": int(odd.degree()),
    }
    if odd.degree() != 4:
        record["status"] = "ODD_PART_NOT_QUARTIC"
        return record

    e, d, c, b, a = [KU(odd[index]) for index in range(5)]
    invariant_i = 12 * a * e - 3 * b * d + c**2
    invariant_j = (
        72 * a * c * e + 9 * b * c * d - 27 * a * d**2
        - 27 * b**2 * e - 2 * c**3
    )
    A4 = -27 * invariant_i
    A6 = -27 * invariant_j
    record["jacobian"] = jacobian_profile(A4, A6)
    record["status"] = "PASS_EXACT_BINARY_QUARTIC_JACOBIAN"
    return record


trials = []
for rho89_sign in (1,):
    for rho68_sign in (1,):
        for chord_sign in (1,):
            for function_sign in (1,):
                rows89, diag89 = toric_rows(
                    F(89), rho89_sign * roots[F(89)]["rho"], chord_sign, function_sign
                )
                rows68, diag68 = toric_rows(
                    F(68), rho68_sign * roots[F(68)]["rho"], chord_sign, function_sign
                )
                matrix_rows, kernel = canonical_kernel(rows89 + rows68)
                compiled = compile_kernel(kernel) if kernel.nrows() == 2 else {
                    "status": "UNEXPECTED_KERNEL_DIMENSION"
                }
                trial = {
                    "rho89_sign": int(rho89_sign),
                    "rho68_sign": int(rho68_sign),
                    "chord_sign": int(chord_sign),
                    "function_sign": int(function_sign),
                    "condition_rows": int(matrix_rows.nrows()),
                    "condition_rank_over_F103": int(matrix_rows.rank()),
                    "kernel_dimension": int(kernel.nrows()),
                    "root89_components": diag89,
                    "root68_components": diag68,
                    "compiled": compiled,
                }
                trials.append(trial)
                profile = compiled.get("jacobian", {})
                print(
                    "H92A5A5Q10TORIC|rho89_sign={}|rho68_sign={}|chord_sign={}|"
                    "function_sign={}|rank={}|kernel={}|odd_degree={}|ADE={}|"
                    "root_rank={}|euler={}|status={}".format(
                        rho89_sign, rho68_sign, chord_sign, function_sign,
                        matrix_rows.rank(), kernel.nrows(), compiled.get("odd_degree"),
                        profile.get("ADE"), profile.get("root_rank"),
                        profile.get("euler_number"), compiled["status"],
                    ),
                    flush=True,
                )

target_trials = [
    trial for trial in trials
    if trial["compiled"].get("jacobian", {}).get("ADE") == ["I0*", "I2", "I4", "I4"]
    and trial["compiled"]["jacobian"]["root_rank"] == 11
    and trial["compiled"]["jacobian"]["euler_number"] == 24
]

inputs = (POOL, SHELL)
payload = {
    "schema": "elkies-k3.q24-2a5-direct-q10-toric-rr-mod103.v1",
    "status": (
        "PASS_EXACT_MOD103_TORIC_Q10_RR_PLANE_WITH_ROOT_JUMP"
        if target_trials else "MOD103_TORIC_Q10_RR_MISMATCH"
    ),
    "prime": int(p),
    "method": {
        "ambient_dimension": 15,
        "smooth_collision_rank": 10,
        "post_collision_dimension": 5,
        "toric_condition_rows": 13,
        "expected_toric_rank": 3,
        "expected_kernel_dimension": 2,
        "enumerated_orientation_sign_trials": len(trials),
        "groebner_basis_used": False,
    },
    "horizontal_section": {
        "MW_word": "2*P_aff-2*P1+3*P1229+2*P230+2*P14",
        "projective_degrees_mod103": [int(X.degree()), int(Y.degree()), int(Z.degree())],
        "Z_coefficients_low_to_high": [int(value) for value in Z.list()],
    },
    "collision_basis": [
        {
            "AA_coefficients_low_to_high": [int(value) for value in AA.list()],
            "BB_coefficients_low_to_high": [int(value) for value in BB.list()],
        }
        for AA, BB in pairs
    ],
    "trials": trials,
    "target_trial_indices": [trials.index(trial) for trial in target_trials],
    "inputs": {
        str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in inputs
    },
    "proof_boundary": (
        "This is an exact finite-field resolved RR and Jacobian-fibre classifier. "
        "The I0*+I2+2I4 specialized profile has root rank 11 and is recorded as "
        "a characteristic-103 root jump relative to the exact generic lattice "
        "target 3A3 of rank 9. It does not prove the characteristic-zero quartic "
        "or fibre profile."
    ),
}
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92A5A5Q10TORIC|trials={}|targets={}|target_indices={}|status={}".format(
        len(trials), len(target_trials), payload["target_trial_indices"], payload["status"]
    ),
    flush=True,
)
if not target_trials:
    raise ArithmeticError("pinned exact mod-103 toric RR profile mismatch")
