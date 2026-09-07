#!/usr/bin/env sage -python
"""Reconstruct the low-pole A11 bridge M from sampled exact trace inputs.

Reconstruct both quartic branches of E5, S7, S14 and S17 individually and
select them by their exact lattice-predicted pole/component profiles.  Then
apply the exact lattice word

    M-a = -AJ(E5_actual) + AJ(S0) + AJ(S1) - 4*Qminus.

The invariant quartic zero is an old spinor component, not the pinned orbit64
zero.  Thus lattice M is represented by M-a, with pole order eight and I12
depth one.  The independently reconstructed pinned zero supplies the stronger
gate O_pinned-M=-M, of pole order five and depth three.  Newton interpolation,
extended-Euclid Pade reconstruction and exact F_p(T) group law are used; no
multivariate solve or Groebner basis is used.
"""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import (
    EllipticCurve,
    GF,
    PolynomialRing,
    PowerSeriesRing,
    QQ,
    ZZ,
    block_diagonal_matrix,
    identity_matrix,
    matrix,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--traces",
    type=Path,
    default=LOCAL / "q24-a11-corrected-bridge-traces-220-mod100003.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-a11-bridge-m-section-marked-mod100003.json",
)
args = parser.parse_args()

TRACES = args.traces.resolve()
A11 = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
ROUTE = GENERATED / "elkies-k3-h3-q24-a11-degree13-e5-bridge-route.json"
PHYSICAL = LOCAL / "q24-d12-orbit42-i8star-physical-marking-qq.json"
NEIGHBOURS = LOCAL / "q24-downstream-lift/d12-c10a-zero-q6-all.json"
PINNED_ZERO = LOCAL / "q24-a11-pinned-zero-section-mod100003.json"
POINTED = LOCAL / "q24-a11-pointed-opposite-section-qq.json"
INPUTS = (TRACES, A11, ROUTE, PHYSICAL, NEIGHBOURS, PINNED_ZERO, POINTED)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

traces = json.loads(TRACES.read_text())
a11 = json.loads(A11.read_text())
route = json.loads(ROUTE.read_text())
physical = json.loads(PHYSICAL.read_text())
neighbours = json.loads(NEIGHBOURS.read_text())
pinned_zero = json.loads(PINNED_ZERO.read_text())
pointed = json.loads(POINTED.read_text())
assert traces["status"] == "PASS_Q24_A11_CORRECTED_BRIDGE_TRACE_SAMPLES_MODP"
assert a11["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert route["status"] == "PASS_EXACT_A11_DEGREE13_E5_BRIDGE_ROUTE"
assert pinned_zero["status"] == "PASS_Q24_A11_PINNED_ZERO_SECTION_RECONSTRUCTION_MODP"
assert pointed["status"] == "PASS_EXACT_A11_POINTED_OPPOSITE_SECTION_QQ"

# The selected physical orientation sends abstract D12 component 6 to C10.
# Its effective negative root is the equation's pointed-quartic zero section
# relative to the pinned orbit64 frame.
selected_physical = next(
    row for row in physical["orientation_candidates"]
    if row["section_meets_physical_components"] == ["C10"]
)
anchor_root_index = selected_physical["abstract_to_physical"].index("C10")
assert anchor_root_index == 6
neighbor = next(row for row in neighbours["neighbors"] if int(row["orbit_index"]) == 64)
transition = block_diagonal_matrix(
    identity_matrix(ZZ, 2), matrix(ZZ, neighbor["child_root_adapted_basis"])
) * matrix(ZZ, neighbor["neighbor_basis"])
anchor_parent = vector(ZZ, [0, 0] + [(-1 if index == anchor_root_index else 0) for index in range(17)])
anchor_child = anchor_parent * transition.inverse().change_ring(ZZ)
anchor_mw = vector(ZZ, anchor_child[-6:])
assert anchor_child[1] == 1 and anchor_mw == vector(ZZ, (1, -2, 1, -1, 0, 0))

p = ZZ(traces["prime"])
F = GF(p)
R = PolynomialRing(F, "T")
T = R.gen()
K = R.fraction_field()


def red(value):
    value = QQ(value)
    if value.denominator() % p == 0:
        raise ZeroDivisionError(f"bad denominator modulo {p}")
    return F(value.numerator()) / F(value.denominator())


A = R([red(value) for value in a11["child"]["minimal_A_coefficients_low_to_high"]])
B = R([red(value) for value in a11["child"]["minimal_B_coefficients_low_to_high"]])

curve_names = ("E5", "S0", "S1", "S7")


def interpolation_polynomial(values):
    interpolation = R.zero()
    modulus = R.one()
    for parameter, value in values:
        scale = modulus(parameter)
        if not scale:
            raise ArithmeticError("duplicate interpolation parameter")
        interpolation += ((value - interpolation(parameter)) / scale) * modulus
        modulus *= T - parameter
    interpolation %= modulus
    return interpolation, modulus


def pade_sequence(values):
    interpolation, modulus = interpolation_polynomial(values)
    r0, r1 = modulus, interpolation
    t0, t1 = R.zero(), R.one()
    candidates = []
    while r1:
        numerator, denominator = r1, t1
        common = numerator.gcd(denominator)
        numerator //= common
        denominator //= common
        if denominator:
            scale = denominator.leading_coefficient()
            numerator /= scale
            denominator /= scale
            if all(
                denominator(parameter)
                and numerator(parameter) == value * denominator(parameter)
                for parameter, value in values
            ):
                candidates.append((numerator, denominator))
        quotient, r2 = r0.quo_rem(r1)
        t2 = t0 - quotient * t1
        r0, r1 = r1, r2
        t0, t1 = t1, t2
    return candidates


def section_pole_order(section):
    X, Y, Z = section
    finite = int(Z.degree())
    x_excess = int(X.degree()) - 2 * finite - 4
    y_excess = int(Y.degree()) - 3 * finite - 6
    infinity = max(0, -((-x_excess) // 2), -((-y_excess) // 3))
    return finite + infinity


def reconstruct_section(samples, pole_order):
    x_candidates = []
    for numerator, denominator in pade_sequence([(tau, x) for tau, x, unused in samples]):
        if denominator.degree() > 2 * pole_order or not denominator.is_square():
            continue
        Z = denominator.sqrt().monic()
        if numerator.degree() > 2 * pole_order + 4:
            continue
        x_candidates.append((numerator, Z))
    answers = []
    y_pade = pade_sequence([(tau, y) for tau, unused, y in samples])
    for X, Z in x_candidates:
        for numerator, denominator in y_pade:
            scale = denominator.leading_coefficient()
            numerator /= scale
            denominator /= scale
            if denominator != Z**3 or numerator.degree() > 3 * pole_order + 6:
                continue
            Y = numerator
            section = (X, Y, Z)
            if (
                Y**2 == X**3 + A * X * Z**4 + B * Z**6
                and section_pole_order(section) == pole_order
            ):
                answers.append((X, Y, Z))
    unique = {}
    for X, Y, Z in answers:
        unique[(tuple(X.list()), tuple(Y.list()), tuple(Z.list()))] = (X, Y, Z)
    return list(unique.values())


def denominator_power_root(value, exponent):
    denominator = R(K(value).denominator()).monic()
    answer = R.one()
    for factor, multiplicity in denominator.factor():
        if int(multiplicity) % exponent:
            return None
        answer *= factor.monic() ** (int(multiplicity) // exponent)
    return answer.monic()


def normalized_point(point):
    if point.is_zero():
        return None
    x, y = map(K, point.xy())
    Zx = denominator_power_root(x, 2)
    Zy = denominator_power_root(y, 3)
    if Zx is None or Zy is None or Zx != Zy:
        return None
    Z = Zx
    X = R(x * Z**2)
    Y = R(y * Z**3)
    if Y**2 != X**3 + A * X * Z**4 + B * Z**6:
        return None
    return X, Y, Z


Delta = R(-16 * (4 * A**3 + 27 * B**2))
i12 = [(factor.monic(), int(exponent)) for factor, exponent in Delta.factor() if int(exponent) == 12]
if len(i12) != 1 or i12[0][0].degree() != 1:
    raise ArithmeticError("good reduction does not retain one rational I12 fibre")
beta = -i12[0][0][0]
PS = PowerSeriesRing(F, "s", default_prec=15)
s = PS.gen()
A_series = PS(A(s + beta))
B_series = PS(B(s + beta))
center = PS(-3 * B_series[0] / (2 * A_series[0]))
for unused in range(7):
    center = (center + (-A_series / 3) / center) / 2
if (center**2 + A_series / 3).valuation() < 14:
    raise ArithmeticError("I12 formal center did not converge")


def component_depth(section):
    X, Y, Z = section
    Zs = PS(Z(s + beta))
    if not Zs:
        return None
    xs = PS(X(s + beta)) / Zs**2
    ys = PS(Y(s + beta)) / Zs**3
    return min(int((xs - center).valuation()), int(ys.valuation()), 6)


def section_row(section):
    X, Y, Z = section
    return {
        "X_coefficients_low_to_high": [int(value) for value in X.list()],
        "Y_coefficients_low_to_high": [int(value) for value in Y.list()],
        "Z_coefficients_low_to_high": [int(value) for value in Z.list()],
        "degrees_X_Y_Z": [int(X.degree()), int(Y.degree()), int(Z.degree())],
        "finite_denominator_degree": int(Z.degree()),
        "P_dot_equation_zero": section_pole_order(section),
        "I12_component_depth_up_to_negation": component_depth(section),
        "exact_modp_weierstrass_identity": True,
    }


def trace_samples(name, branch):
    values = []
    for sample in traces["samples"]:
        tau = F(sample["tau"])
        row = next(
            item
            for item in sample["points"]
            if item["curve"] == name and int(item["quartic_branch_sign"]) == branch
        )
        values.append((tau, F(row["AJ_x"]), F(row["AJ_y"])))
    return values


expected_profiles = {
    "E5": (30, 3),
    "S0": (6, 1),
    "S1": (3, 3),
    "S7": (0, 3),
}
trace_sections = {}
orientation_audit = []
for name in curve_names:
    expected_pole, expected_depth = expected_profiles[name]
    selected = []
    for branch in (1, -1):
        samples = trace_samples(name, branch)
        answers = reconstruct_section(samples, expected_pole)
        accepted = [section for section in answers if component_depth(section) == expected_depth]
        orientation_audit.append(
            {
                "curve": name,
                "quartic_branch_sign": branch,
                "usable_samples": len(samples),
                "expected_profile": [expected_pole, expected_depth],
                "reconstruction_count": len(answers),
                "profile_compatible_count": len(accepted),
            }
        )
        selected.extend((branch, section) for section in accepted)
    if len(selected) != 1:
        raise ArithmeticError(
            f"{name}: expected one profile-compatible trace branch, got {len(selected)}"
        )
    trace_sections[name] = selected[0]

generic_curve = EllipticCurve(K, [0, 0, 0, K(A), K(B)])
trace_points = {
    name: generic_curve(
        K(section[0]) / K(section[2] ** 2),
        K(section[1]) / K(section[2] ** 3),
    )
    for name, (unused_branch, section) in trace_sections.items()
}


def reduced_section(row):
    return tuple(
        R([red(value) for value in row[key]])
        for key in (
            "X_coefficients_low_to_high",
            "Y_coefficients_low_to_high",
            "Z_coefficients_low_to_high",
        )
    )


Q_section = reduced_section(pointed["section"])
Q_generic = generic_curve(
    K(Q_section[0]) / K(Q_section[2] ** 2),
    K(Q_section[1]) / K(Q_section[2] ** 3),
)
zero_section = reduced_section(pinned_zero["section"])
zero_generic = generic_curve(
    K(zero_section[0]) / K(zero_section[2] ** 2),
    K(zero_section[1]) / K(zero_section[2] ** 3),
)

# The profile of S7 is invariant under inversion.  Its pairing with the
# already pinned Qminus fixes the actual MW sign used by the final gate.
s7_point = trace_points["S7"]
q_plus_s7 = normalized_point(Q_generic + s7_point)
q_minus_s7 = normalized_point(Q_generic - s7_point)
q_s7_profiles = {
    "QplusS7": [section_pole_order(q_plus_s7), component_depth(q_plus_s7)],
    "QminusS7": [section_pole_order(q_minus_s7), component_depth(q_minus_s7)],
}
if q_s7_profiles == {"QplusS7": [0, 1], "QminusS7": [1, 5]}:
    s7_trace_sign = 1
elif q_s7_profiles == {"QplusS7": [1, 5], "QminusS7": [0, 1]}:
    s7_trace_sign = -1
    s7_point = -s7_point
else:
    raise ArithmeticError(f"could not pin S7 against Qminus: {q_s7_profiles}")

construction_names = ("E5", "S0", "S1")
sign_audit = []
marked_solutions = []
for signs in itertools.product((1, -1), repeat=len(construction_names)):
    signed = dict(zip(construction_names, signs))
    candidate = (
        -signed["E5"] * trace_points["E5"]
        + signed["S0"] * trace_points["S0"]
        + signed["S1"] * trace_points["S1"]
        - 4 * Q_generic
    )
    section = normalized_point(candidate)
    section_profile = None if section is None else [
        section_pole_order(section), component_depth(section)
    ]
    difference = None
    difference_profile = None
    reference_profiles = None
    if section_profile == [8, 1]:
        difference = normalized_point(zero_generic - candidate)
        if difference is not None:
            difference_profile = [
                section_pole_order(difference), component_depth(difference)
            ]
        reference = normalized_point(candidate - s7_point)
        reference_profiles = {
            "bridge_minus_S7": None if reference is None else [
                section_pole_order(reference), component_depth(reference)
            ]
        }
    sign_audit.append(
        {
            "trace_signs": signed,
            "bridge_profile": section_profile,
            "Opinned_minus_bridge_profile": difference_profile,
            "reference_profiles": reference_profiles,
        }
    )
    if (
        section_profile == [8, 1]
        and difference_profile == [5, 3]
        and reference_profiles["bridge_minus_S7"] == [11, 4]
    ):
        marked_solutions.append((signed, candidate, section, difference))
if len(marked_solutions) != 1:
    raise ArithmeticError(
        f"expected one triply marked corrected bridge, got {len(marked_solutions)}; "
        f"audit={sign_audit}"
    )
selected_signs, bridge_generic, bridge_section, zero_minus_bridge = marked_solutions[0]
X, Y, Z = bridge_section
orientation = {
    name: ("branch_plus" if trace_sections[name][0] == 1 else "branch_minus")
    for name in curve_names
}
orientation["trace_signs"] = selected_signs
orientation["S7_marking_trace_sign"] = s7_trace_sign
orientation["Q_coefficient_on_stored_positive_section"] = -4

payload = {
    "schema": "elkies-k3.h3-q24-a11-bridge-m-section-marked-modp.v2",
    "status": "PASS_Q24_A11_BRIDGE_M_SECTION_MARKED_RECONSTRUCTION_MODP",
    "prime": int(p),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "superseded_bridge_word": route["bridge_word"],
    "corrected_equation_word": {
        "formula": "M-a=-AJ(E5_actual)+AJ(S0)+AJ(S1)-4*Qminus",
        "actual_E5_equation_MW": [0, 6, -1, -1, 0, -1],
        "target_equation_MW": [0, 2, -1, 1, 0, 1],
        "verified_modp_group_law": True,
    },
    "selected_orientation": orientation,
    "equation_zero_translation": {
        "physical_component": "C10",
        "abstract_parent_root_index": anchor_root_index,
        "pinned_A11_MW": [int(value) for value in anchor_mw],
        "relation": "equation_point_for_pinned_M_has_MW=M-a",
        "profile_of_M_minus_equation_zero": [8, 1],
    },
    "orientation_audit": orientation_audit,
    "sign_audit": sign_audit,
    "section": {
        **section_row(bridge_section),
        "pinned_lattice_P_dot_O": 5,
        "A11_MW_Abel_Jacobi": [1, 0, 0, 0, 0, 1],
    },
    "independent_marking_gate": {
        "relation": "O_pinned(eq)-M(eq)=-M",
        "section": section_row(zero_minus_bridge),
        "expected_profile": [5, 3],
        "second_reference": {
            "relation": "M(eq)-S7(eq)",
            "expected_profile": [11, 4],
        },
    },
    "method": "individual univariate Newton/Pade trace reconstruction followed by exact F_p(T) group law and independent pinned-zero marking",
    "large_Groebner_required": False,
    "proof_boundary": (
        "Exact over the pinned good-reduction field with an exact characteristic-zero lattice word. "
        "The displayed equation bridge still requires characteristic-zero coefficient lifting and literal QQ verification."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11BRIDGEMMODP|prime={}|orientation={}|degrees={},{},{}|equation_PO={}|lattice_PO=5|status={}".format(
        p,
        ",".join(f"{name}:{orientation[name]}" for name in curve_names)
        + f",Qcoef:{orientation['Q_coefficient_on_stored_positive_section']}",
        X.degree(),
        Y.degree(),
        section_pole_order(bridge_section),
        Z.degree(),
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
