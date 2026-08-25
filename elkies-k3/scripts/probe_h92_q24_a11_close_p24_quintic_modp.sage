#!/usr/bin/env sage -python
"""Probe ``close_P24`` in the selected R3-zero A11 marking modulo a good prime.

This composes the two already-certified binary-quartic neighbour maps:

    q8/D13 --q24--> D12 --q6--> A11.

The old exact q24 section lies in the base locus of the first resolved RR
pencil.  This script deliberately tests the tempting naive restriction via
the two pointed-quartic signs before the second pencil.  Exact NS replay in
the compatible R3-zero marking gives degree 46; any other degrees prove that
the strict transform must first be evaluated in the resolved basepoint chart.
This is a bounded modular implementation audit, not a characteristic-zero
section or trace certificate.  It uses only univariate arithmetic and exact
elliptic/quartic formulas, never a Groebner basis.
"""

import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, is_prime, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--trace-tau", type=int, default=2)
parser.add_argument("--trace-samples", type=int, default=1)
parser.add_argument(
    "--trace-branch",
    choices=("pole_plus", "pole_minus", "both"),
    default="both",
)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-a11-close-p24-quintic-modp.json",
)
args = parser.parse_args()

p = ZZ(args.prime)
if not is_prime(p) or p in (2, 3):
    raise SystemExit("--prime must be a prime other than 2 or 3")

PARENT_PATH = LOCAL / "q8-corrected2cover-qq-child.json"
SECTION_PATH = LOCAL / "q8-q24-horizontal-section-qq.json"
Q24_PATH = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
A11_PATH = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
CANDIDATES_PATH = LOCAL / "q24-orbit42-exact-section-candidates-qq.json"
MODEL_PATH = LOCAL / "q24-orbit42-zero-pole-seeds-mod-43.json"
INPUTS = (
    PARENT_PATH,
    SECTION_PATH,
    Q24_PATH,
    A11_PATH,
    CANDIDATES_PATH,
    MODEL_PATH,
)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

parent = json.loads(PARENT_PATH.read_text())
section = json.loads(SECTION_PATH.read_text())
q24 = json.loads(Q24_PATH.read_text())
a11 = json.loads(A11_PATH.read_text())
candidates = json.loads(CANDIDATES_PATH.read_text())
model = json.loads(MODEL_PATH.read_text())["exact_model"]

assert parent["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"
assert section["status"] == "PASS_EXACT_Q24_HORIZONTAL_SECTION"
assert q24["status"] == "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR"
assert a11["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert candidates["status"] == "PASS_EXACT_Q42_ORBIT42_SECTION_CANDIDATES_QQ"

F = GF(p)
UR = PolynomialRing(F, "U")
U = UR.gen()
KU = UR.fraction_field()
VR = PolynomialRing(F, "V")
V = VR.gen()
KV = VR.fraction_field()


def red(value):
    value = QQ(value)
    if value.denominator() % p == 0:
        raise ZeroDivisionError(f"bad denominator modulo {p}")
    return F(value.numerator()) / F(value.denominator())


def red_poly(values, ring):
    return ring([red(value) for value in values])


def reduce_qq_rational(text, name):
    """Parse a QQ(name) string and reduce it into the matching F_p field."""

    source_ring = PolynomialRing(QQ, name)
    source_field = source_ring.fraction_field()
    value = source_field(str(text))
    target_ring = UR if name == "U" else VR
    target_field = target_ring.fraction_field()
    numerator = target_ring([red(value.numerator()[i]) for i in range(value.numerator().degree() + 1)])
    denominator = target_ring([red(value.denominator()[i]) for i in range(value.denominator().degree() + 1)])
    return target_field(numerator) / target_field(denominator)


def map_degree(value, ring):
    value = ring.fraction_field()(value)
    return max(
        int(ring(value.numerator()).degree()),
        int(ring(value.denominator()).degree()),
    )


def eval_v(value, argument):
    value = KV(value)
    return KU(UR(value.numerator())(argument)) / KU(UR(value.denominator())(argument))


def eval_u(value, argument):
    value = KU(value)
    return KV(VR(value.numerator())(argument)) / KV(VR(value.denominator())(argument))


def square_root_rational(value, ring):
    value = ring.fraction_field()(value)
    if not value.is_square():
        raise ArithmeticError("expected rational-function square")
    answer = value.sqrt()
    assert answer**2 == value
    return answer


def covariants_at(coefficients, x_value):
    e, d, c, b, a = coefficients
    g0 = b**2 / 16 - a * c / 6
    g1 = b * c / 12 - a * d / 2
    g2 = c**2 / 12 - b * d / 8 - a * e
    g3 = c * d / 12 - b * e / 2
    g4 = d**2 / 16 - c * e / 6
    ux = 4 * a * x_value**3 + 3 * b * x_value**2 + 2 * c * x_value + d
    uy = b * x_value**3 + 2 * c * x_value**2 + 3 * d * x_value + 4 * e
    g = g0 * x_value**4 + g1 * x_value**3 + g2 * x_value**2 + g3 * x_value + g4
    gx = 4 * g0 * x_value**3 + 3 * g1 * x_value**2 + 2 * g2 * x_value + g3
    gy = g1 * x_value**3 + 2 * g2 * x_value**2 + 3 * g3 * x_value + 4 * g4
    h = (ux * gy - uy * gx) / 8
    return g, h


# -------------------------------------------------------------------------
# 1. Restrict the exact q24 pencil to its old rational section.
# -------------------------------------------------------------------------
old_child = parent["child"]
A13 = red_poly(old_child["minimal_A_coefficients_low_to_high"], UR)
B13 = red_poly(old_child["minimal_B_coefficients_low_to_high"], UR)
sec = section["section"]
Z24 = red_poly(sec["Z_coefficients_low_to_high"], UR)
X24 = red_poly(sec["X_coefficients_low_to_high"], UR)
Y24 = red_poly(sec["Y_coefficients_low_to_high"], UR)
assert Y24**2 == X24**3 + A13 * X24 * Z24**4 + B13 * Z24**6

plane = q24["rr"]["plane_2x56"]
assert len(plane) == 2 and all(len(row) == 56 for row in plane)
pairs24 = []
for row in plane:
    AA = red_poly(row[:41], UR)
    BB = red_poly(row[41:], UR)
    pairs24.append((KU(AA) / KU(Z24**2), KU(BB) / KU(Z24)))

(a240, b240), (a241, b241) = pairs24
# At the old section itself the chord slope is infinite, so the restricted
# pencil value is the ratio of the two slope coefficients.
V_on_curve = KU(b241 / b240)
x_section = KU(X24) / KU(Z24**2)
y_section = KU(Y24) / KU(Z24**3)
# The chord discriminant was compiled with marked point ``(x,-y)``.  Test
# both tangent signs explicitly: this is a convention diagnostic, not an
# assumption about which branch represents the strict transform through the
# resolved RR base locus.
tangent_slopes = {
    "opposite_tangent": -(3 * x_section**2 + KU(A13)) / (2 * y_section),
    "section_tangent": (3 * x_section**2 + KU(A13)) / (2 * y_section),
}
tangent_maps = {
    name: KU((a241 + b241 * slope) / (a240 + b240 * slope))
    for name, slope in tangent_slopes.items()
}

q24_coefficients = [
    reduce_qq_rational(text, "V")
    for text in q24["quartic"]["coefficients_in_U_low_to_high"]
]
q24_value = sum(eval_v(coefficient, V_on_curve) * U**i for i, coefficient in enumerate(q24_coefficients))
w24 = square_root_rational(q24_value, UR)
q24_tangent_values = {
    name: KU(sum(
        eval_v(coefficient, Vmap) * U**i
        for i, coefficient in enumerate(q24_coefficients)
    ))
    for name, Vmap in tangent_maps.items()
}


def chord_discriminant(x_point, y_point, old_a, slope):
    return (
        slope**4 - 6 * x_point * slope**2 + 8 * y_point * slope
        - 3 * x_point**2 - 4 * old_a
    )


direct_tangent_discriminants = {
    name: KU(chord_discriminant(x_section, -y_section, KU(A13), slope))
    for name, slope in tangent_slopes.items()
}
recovered_tangent_slopes = {
    name: KU((a241 - Vmap * a240) / (Vmap * b240 - b241))
    for name, Vmap in tangent_maps.items()
}
assert all(recovered_tangent_slopes[name] == tangent_slopes[name] for name in tangent_maps)

A12 = red_poly(q24["child"]["minimal_A_coefficients_low_to_high"], VR)
B12 = red_poly(q24["child"]["minimal_B_coefficients_low_to_high"], VR)
raw_A12 = reduce_qq_rational(q24["jacobian_raw"]["A"], "V")
raw_B12 = reduce_qq_rational(q24["jacobian_raw"]["B"], "V")

# Point the q24 quartic at the old I9* branch point.  This is the degree-one
# quartic-to-Jacobian isomorphism; using the invariant covariant here would
# double the fibrewise class and destroy the expected downstream degree.
i9 = next(item for item in old_child["finite_fibres"] if item["kodaira"] == "I9*")
i9_source = PolynomialRing(QQ, "U")(str(i9["factor"]))
alpha24 = red(-i9_source[0] / i9_source[1])
qalpha = sum(q24_coefficients[i] * alpha24**i for i in range(5))
q0 = square_root_rational(qalpha, VR)
qa0, qa1, qa2, qa3, qa4 = q24_coefficients
a24 = qa4
b24 = qa3 + 4 * alpha24 * qa4
c24 = qa2 + 3 * alpha24 * qa3 + 6 * alpha24**2 * qa4
d24 = qa1 + 2 * alpha24 * qa2 + 3 * alpha24**2 * qa3 + 4 * alpha24**3 * qa4

cA = raw_A12 / KV(A12)
cB = raw_B12 / KV(B12)
scale2 = cB / cA
scale3 = square_root_rational(cB, VR)
assert scale2**2 == cA and scale3**2 == cB


def q24_point_to_minimal(Vmap, Wvalue):
    uu = KU(U - alpha24)
    qt = eval_v(q0, Vmap)
    at = eval_v(a24, Vmap)
    bt = eval_v(b24, Vmap)
    ct = eval_v(c24, Vmap)
    dt = eval_v(d24, Vmap)
    xg = (2 * qt * (Wvalue + qt) + dt * uu) / uu**2
    yg = (
        4 * qt**2 * (Wvalue + qt)
        + 2 * qt * (dt * uu + ct * uu**2)
        - dt**2 * uu**2 / (2 * qt)
    ) / uu**3
    a1 = dt / qt
    a2 = ct - dt**2 / (4 * qt**2)
    a3 = 2 * qt * bt
    b2 = a1**2 + 4 * a2
    xraw = 9 * (xg + b2 / 12)
    yraw = 27 * (yg + (a1 * xg + a3) / 2)
    assert yraw**2 == (
        xraw**3
        + eval_v(raw_A12, Vmap) * xraw
        + eval_v(raw_B12, Vmap)
    )
    xmin = KU(xraw / eval_v(scale2, Vmap))
    ymin = KU(yraw / eval_v(scale3, Vmap))
    assert ymin**2 == (
        xmin**3
        + eval_v(KV(A12), Vmap) * xmin
        + eval_v(KV(B12), Vmap)
    )
    return xmin, ymin


candidate_D12_points = [
    ("pole_plus", V_on_curve, *q24_point_to_minimal(V_on_curve, w24)),
    ("pole_minus", V_on_curve, *q24_point_to_minimal(V_on_curve, -w24)),
]
for name, Vmap in tangent_maps.items():
    if q24_tangent_values[name] == 0:
        candidate_D12_points.append(
            (name, Vmap, *q24_point_to_minimal(Vmap, KU.zero()))
        )


def reduce_modulus(value, modulus):
    value = KU(value)
    numerator = UR(value.numerator())
    denominator = UR(value.denominator())
    if denominator.gcd(modulus).degree() != 0:
        raise ZeroDivisionError("multisection denominator meets trace fibre")
    return UR((numerator * denominator.inverse_mod(modulus)) % modulus)


def newton_power_sums(monic):
    degree = monic.degree()
    assert monic[degree] == 1
    sums = [F(degree)]
    for exponent in range(1, degree):
        total = F(exponent) * monic[degree - exponent]
        for index in range(1, exponent):
            total += monic[degree - index] * sums[exponent - index]
        sums.append(-total)
    return sums


def trace_degree14(branch, Vmap, xmap, ymap, tau_integer):
    tau = F(tau_integer)
    equation = UR(Vmap.numerator() - tau * Vmap.denominator())
    common = equation.gcd(equation.derivative())
    if common.degree() != 0:
        raise ArithmeticError(f"{branch}: trace fibre is not etale")
    modulus = equation.monic()
    if modulus.degree() != 14:
        raise ArithmeticError(f"{branch}: trace fibre degree {modulus.degree()}, expected 14")
    x_residue = reduce_modulus(xmap, modulus)
    y_residue = reduce_modulus(ymap, modulus)
    A_value = A12(tau)
    B_value = B12(tau)
    if (y_residue**2 - x_residue**3 - A_value * x_residue - B_value) % modulus:
        raise ArithmeticError(f"{branch}: D12 point identity failed modulo trace fibre")

    powers = [UR.one()]
    for unused in range(7):
        powers.append((powers[-1] * x_residue) % modulus)
    columns = powers + [(y_residue * powers[index]) % modulus for index in range(7)]
    assert len(columns) == 15
    evaluation = matrix(F, 14, 15, lambda row, column: columns[column][row])
    kernel = evaluation.right_kernel().basis_matrix()
    if kernel.nrows() != 1:
        raise ArithmeticError(f"{branch}: L(15O) trace kernel dimension {kernel.nrows()}")
    relation = kernel[0]

    XR = PolynomialRing(F, "Xtrace")
    Xtrace = XR.gen()
    afun = sum(relation[index] * Xtrace**index for index in range(8))
    bfun = sum(relation[8 + index] * Xtrace**index for index in range(7))
    intersection = afun**2 - (Xtrace**3 + A_value * Xtrace + B_value) * bfun**2
    if intersection.degree() != 15:
        raise ArithmeticError(f"{branch}: residual intersection degree {intersection.degree()}")
    residual_root_sum = -intersection[14] / intersection[15]
    power_sums = newton_power_sums(modulus)
    trace_x = sum(x_residue[index] * power_sums[index] for index in range(14))
    residual_x = residual_root_sum - trace_x
    residual_b = bfun(residual_x)
    if not residual_b:
        raise ArithmeticError(f"{branch}: residual B value vanished")
    residual_y = -afun(residual_x) / residual_b
    curve = EllipticCurve(F, [0, 0, 0, A_value, B_value])
    aj = -curve(residual_x, residual_y)
    return {
        "branch": branch,
        "tau": int(tau),
        "degree": 14,
        "trace_space": "L(15O)",
        "AJ_x": int(aj[0]),
        "AJ_y": int(aj[1]),
    }


selected_trace_branches = {
    "both": {"pole_plus", "pole_minus"},
    "pole_plus": {"pole_plus"},
    "pole_minus": {"pole_minus"},
}[args.trace_branch]
trace_maps = [
    (branch, Vmap, x12, y12)
    for branch, Vmap, x12, y12 in candidate_D12_points
    if branch in selected_trace_branches
]
degree14_traces = []
tau_integer = args.trace_tau
while len(degree14_traces) < args.trace_samples * len(trace_maps):
    if tau_integer >= p:
        raise ArithmeticError("ran out of finite-field trace specializations")
    batch = []
    try:
        for branch, Vmap, x12, y12 in trace_maps:
            batch.append(trace_degree14(branch, Vmap, x12, y12, tau_integer))
    except (ArithmeticError, ZeroDivisionError, ValueError):
        tau_integer += 1
        continue
    degree14_traces.extend(batch)
    tau_integer += 1
print(
    "Q24CLOSEP24_AJ14|prime={}|samples={}|branches={}|first={}|status=PASS".format(
        p,
        args.trace_samples,
        args.trace_branch,
        ";".join(
            f"{row['branch']}@{row['tau']}:{row['AJ_x']},{row['AJ_y']}"
            for row in degree14_traces[: len(trace_maps)]
        ),
    ),
    flush=True,
)

print(
    "A11CLOSEP24_Q24|"
    f"prime={p}|pole_V_degree={map_degree(V_on_curve, UR)}|"
    "tangent_V_degrees="
    f"{','.join(name + ':' + str(map_degree(Vmap, UR)) for name, Vmap in tangent_maps.items())}|"
    "tangent_direct_discriminants="
    f"{','.join(name + ':' + str(int(value == 0)) for name, value in direct_tangent_discriminants.items())}|"
    "tangent_quartic_zeros="
    f"{','.join(name + ':' + str(int(q24_tangent_values[name] == 0)) for name in tangent_maps)}|"
    f"quartic_square=1|pointed_D12_points={len(candidate_D12_points)}|status=PASS",
    flush=True,
)

# -------------------------------------------------------------------------
# 2. Reconstruct the exact orbit42 pencil on the minimal D12 model.
# -------------------------------------------------------------------------
u_of_V = reduce_qq_rational(a11["coordinate_change"]["u_of_V"], "V")
x_scale = reduce_qq_rational(a11["coordinate_change"]["x_scale"], "V")
y_scale = reduce_qq_rational(a11["coordinate_change"]["y_scale"], "V")

selected = candidates["candidates"][a11["selected_section"]["candidate_index"]]
uX = red_poly(selected["X_coefficients_low_to_high"], UR)
uY = red_poly(selected["Y_coefficients_low_to_high"], UR)
uZ = red_poly(selected["Z_coefficients_low_to_high"], UR)
x_parent = x_scale * eval_u(KU(uX) / KU(uZ**2), u_of_V)
y_parent = y_scale * eval_u(KU(uY) / KU(uZ**3), u_of_V)
assert y_parent**2 == x_parent**3 + KV(A12) * x_parent + KV(B12)

denominator = VR(x_parent.denominator()).monic()
Z_parent = VR.one()
for factor, multiplicity in denominator.factor():
    if int(multiplicity) % 2:
        raise ArithmeticError("orbit42 parent denominator is not a square")
    Z_parent *= factor.monic() ** (int(multiplicity) // 2)
X_parent = VR(x_parent * Z_parent**2)
Y_parent = VR(y_parent * Z_parent**3)

alpha = red(model["I8star_root"])
modulus = Z_parent**2
X_inverse = VR(X_parent.inverse_mod(modulus))
pairs11 = []
for BB in (VR.one(), V):
    AA = VR((BB * Y_parent * X_inverse) % modulus)
    AA -= AA(alpha) / Z_parent(alpha) ** 2 * Z_parent**2
    pairs11.append((KV(AA) / KV(Z_parent**2), KV(BB) / KV(Z_parent)))

(a110, b110), (a111, b111) = pairs11
base_candidates = []
for branch, Vmap, x12, y12 in candidate_D12_points:
    x_parent_on_curve = eval_v(x_parent, Vmap)
    y_parent_on_curve = eval_v(y_parent, Vmap)
    chord = (y12 + y_parent_on_curve) / (x12 - x_parent_on_curve)
    new_base = (
        eval_v(a111, Vmap) + eval_v(b111, Vmap) * chord
    ) / (
        eval_v(a110, Vmap) + eval_v(b110, Vmap) * chord
    )
    new_base = KU(new_base)
    base_candidates.append((map_degree(new_base, UR), branch, Vmap, new_base))

expected_degree = 46
resolved_degree = [row for row in base_candidates if row[0] == expected_degree]
if len(resolved_degree) != 1:
    payload = {
        "schema": "elkies-k3.h3-q24-a11-close-p24-marking-audit-modp.v2",
        "status": "Q24_A11_CLOSE_P24_NEEDS_RESOLVED_BASEPOINT_RESTRICTION",
        "prime": int(p),
        "q24_stage": {
            "naive_base_ratio_degree": map_degree(V_on_curve, UR),
            "quartic_square": True,
            "pointed_D12_point_identities": True,
            "tangent_diagnostics": {
                name: {
                    "base_ratio_degree": map_degree(tangent_maps[name], UR),
                    "direct_chord_discriminant_zero": direct_tangent_discriminants[name] == 0,
                    "stored_quartic_zero": q24_tangent_values[name] == 0,
                }
                for name in tangent_maps
            },
            "degree14_AJ_traces": degree14_traces,
        },
        "A11_stage": {
            "naive_q24_pointed_sign_A11_degrees": [row[0] for row in base_candidates],
            "expected_strict_transform_degree_from_NS": expected_degree,
            "shortcut_matches_expected_degree": False,
        },
        "conclusion": (
            "close_P24 meets the resolved q24 RR base locus. Substituting the "
            "infinite-chord coefficient ratio and either pointed-quartic sign "
            "does not represent its strict transform. Evaluate the curve along "
            "the certified blow-up chart before applying the A11 pencil. The "
            "compatible R3-zero marking makes this a degree-46 carrier, not the "
            "previously claimed quintic."
        ),
        "proof_boundary": (
            "Exact arithmetic over the declared finite field rejects only the "
            "naive unresolved composition. It does not compute the degree-46 "
            "strict transform or its Abel-Jacobi trace."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "A11CLOSEP24_RESOLUTION_GATE|"
        f"prime={p}|naive_V_degree={map_degree(V_on_curve, UR)}|"
        f"naive_A11_degrees={','.join(str(row[0]) for row in base_candidates)}|"
        f"expected_strict_degree={expected_degree}|"
        f"status={payload['status']}",
        flush=True,
    )
    print(f"OUTPUT|{args.output}", flush=True)
    raise SystemExit(0)
unused_degree, selected_q24_branch, selected_V_map, T_on_curve = resolved_degree[0]

a11_coefficients = []
for text in a11["quartic"]["coefficients_in_T_low_to_high"]:
    source_ring = PolynomialRing(QQ, "T")
    source = source_ring(str(text))
    coefficient = KU.zero()
    for i in range(source.degree() + 1):
        coefficient += KU(red(source[i])) * T_on_curve**i
    a11_coefficients.append(coefficient)

a11_value = sum(
    coefficient * selected_V_map**i
    for i, coefficient in enumerate(a11_coefficients)
)
w11 = square_root_rational(a11_value, UR)

print(
    "A11CLOSEP24_RESOLVED|"
    f"prime={p}|q24_pointed_branch={selected_q24_branch}|"
    f"sign_degrees={','.join(str(row[0]) for row in base_candidates)}|"
    f"A11_degree={expected_degree}|quartic_square=1|"
    f"w_degree={map_degree(w11, UR)}|status=PASS",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h3-q24-a11-close-p24-marking-audit-modp.v2",
    "status": "PASS_A11_CLOSE_P24_DEGREE46_MODP",
    "prime": int(p),
    "q24_stage": {
        "new_base_degree": map_degree(V_on_curve, UR),
        "quartic_square": True,
        "D12_covariant_point_identity": True,
    },
    "A11_stage": {
        "q24_pointed_branch": selected_q24_branch,
        "q24_pointed_branch_A11_degrees": {
            row[1]: row[0] for row in base_candidates
        },
        "new_base_degree": map_degree(T_on_curve, UR),
        "quartic_square": True,
        "quartic_ordinate_degree": map_degree(w11, UR),
    },
    "method": (
        "Composition of the two exact RR pencils and binary-quartic covariants, "
        "reduced modulo the declared good prime; univariate arithmetic only."
    ),
    "proof_boundary": (
        "This verifies the resolved composition and degree-46/square identities "
        "over one finite field. It does not compute the L(6O) trace, orient the "
        "A11 Abel-Jacobi point, lift coordinates to QQ, or prove the next equation."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{args.output}", flush=True)
