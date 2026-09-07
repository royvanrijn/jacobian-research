#!/usr/bin/env sage -python
"""Reconstruct the two A11 q8 horizontal points over one good prime.

The invariant-quartic A11 equation uses C10 as zero.  The q8 divisor is

    O_pinned + P12 - 2F,

so its horizontal points have equation-side MW vectors ``-a`` and ``P12-a``.
Use the exact low-trace words and the sampled D12-to-A11 traces to recover:

* ``O_pinned`` with pole order 3 and I12 depth 2;
* ``P12`` with pole order 11 and I12 depth 4.

The target is additionally required to differ from ``O_pinned`` by a point of
pole order 6 and central I12 depth 6, which is the original pinned q8 section
profile.  Reconstruction is univariate Newton/Padé plus elliptic group law;
no multivariate solve or Groebner basis is used.
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
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--traces",
    type=Path,
    default=LOCAL / "q24-a11-q8-zero-target-traces-mod100003.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-a11-q8-horizontal-points-mod100003.json",
)
parser.add_argument(
    "--zero-output",
    type=Path,
    default=LOCAL / "q24-a11-pinned-zero-section-mod100003.json",
)
parser.add_argument("--stop-after-zero", action="store_true")
parser.add_argument(
    "--bridge-audit-output",
    type=Path,
    default=GENERATED / "elkies-k3-h3-q24-a11-bridge-m-marking-rejection-mod100003.json",
)
args = parser.parse_args()

TRACES = args.traces.resolve()
A11 = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
ROUTE = GENERATED / "elkies-k3-h3-q24-a11-q8-zero-translation-route.json"
BRIDGE_M = LOCAL / "q24-a11-bridge-m-section-marked-qq.json"
POINTED = LOCAL / "q24-a11-pointed-opposite-section-qq.json"
INPUTS = (TRACES, A11, ROUTE, BRIDGE_M, POINTED)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

traces = json.loads(TRACES.read_text())
a11 = json.loads(A11.read_text())
route = json.loads(ROUTE.read_text())
bridge_m = json.loads(BRIDGE_M.read_text())
pointed = json.loads(POINTED.read_text())
assert traces["status"] == "PASS_Q24_A11_Q8_ZERO_TARGET_TRACE_SAMPLES_MODP"
assert a11["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert route["status"] == "PASS_EXACT_A11_Q8_ZERO_TRANSLATION_ROUTE"
assert bridge_m["status"] == "PASS_EXACT_Q24_A11_BRIDGE_M_SECTION_MARKED_QQ"
assert pointed["status"] == "PASS_EXACT_A11_POINTED_OPPOSITE_SECTION_QQ"

p = ZZ(traces["prime"])
F = GF(p)
R = PolynomialRing(F, "T")
T = R.gen()
K = R.fraction_field()
A = R([
    F(QQ(value).numerator()) / F(QQ(value).denominator())
    for value in a11["child"]["minimal_A_coefficients_low_to_high"]
])
B = R([
    F(QQ(value).numerator()) / F(QQ(value).denominator())
    for value in a11["child"]["minimal_B_coefficients_low_to_high"]
])


def interpolation_polynomial(values):
    interpolation = R.zero()
    modulus = R.one()
    for parameter, value in values:
        scale = modulus(parameter)
        if not scale:
            raise ArithmeticError("duplicate interpolation parameter")
        interpolation += ((value - interpolation(parameter)) / scale) * modulus
        modulus *= T - parameter
    return interpolation % modulus, modulus


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


def section_pole_order(section):
    """Total intersection with O, including a possible pole at T=infinity."""
    X, Y, Z = section
    finite = int(Z.degree())
    x_excess = int(X.degree()) - 2 * finite - 4
    y_excess = int(Y.degree()) - 3 * finite - 6
    infinity = max(0, -((-x_excess) // 2), -((-y_excess) // 3))
    return finite + infinity


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


# Formal I12 component-depth evaluator on the exact mod-p section.
Delta = R(-16 * (4 * A**3 + 27 * B**2))
i12 = [(factor.monic(), int(exponent)) for factor, exponent in Delta.factor() if int(exponent) == 12]
if len(i12) != 1 or i12[0][0].degree() != 1:
    raise ArithmeticError("good reduction does not retain one rational I12 fibre")
i12_factor = i12[0][0]
beta = -i12_factor[0]
PS = PowerSeriesRing(F, "s", default_prec=15)
s = PS.gen()
A_series = PS(A(s + beta))
B_series = PS(B(s + beta))
center = PS(-3 * B_series[0] / (2 * A_series[0]))
for unused in range(7):
    center = (center + (-A_series / 3) / center) / 2
if (center**2 + A_series / 3).valuation() < 14:
    raise ArithmeticError("I12 formal center did not converge")
if (center**3 + A_series * center + B_series).valuation() != 12:
    raise ArithmeticError("I12 centered cubic has wrong valuation")


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


def sample_points(sample):
    tau = F(sample["tau"])
    curve = EllipticCurve(F, [0, 0, 0, A(tau), B(tau)])
    points = {
        (row["curve"], int(row["quartic_branch_sign"])): curve(
            F(row["AJ_x"]), F(row["AJ_y"])
        )
        for row in sample["points"]
    }
    return tau, curve, points


# Reconstruct O_pinned from the degree-weight-zero word
# -a=-6*Qminus+7*S5-13*S7+S13-5*S17.  Zero degree weight is essential here:
# it makes the pinned lattice word represent the same class after moving the
# equation zero to C10.  The exact bridge M supplies an independent marking:
# O_pinned(eq)-M(eq)=-M has pole order 5 and I12 depth 3.
generic_curve = EllipticCurve(K, [0, 0, 0, K(A), K(B)])


def reconstruct_individual_trace(name, branch, maximum_pole_order=15):
    samples = []
    for sample in traces["samples"]:
        tau, unused_curve, points = sample_points(sample)
        point = points[(name, branch)]
        if point.is_zero():
            continue
        x, y = point.xy()
        samples.append((tau, F(x), F(y)))
    for pole_order in range(maximum_pole_order + 1):
        answers = reconstruct_section(samples, pole_order)
        if answers:
            if len(answers) != 1:
                raise ArithmeticError(
                    f"{name} branch {branch} has {len(answers)} minimal reconstructions"
                )
            return answers[0]
    raise ArithmeticError(
        f"{name} branch {branch} was not reconstructed through P.O={maximum_pole_order}"
    )


trace_sections = {
    (name, branch): reconstruct_individual_trace(name, branch)
    for name in ("S5", "S7", "S13", "S17")
    for branch in (1, -1)
}
trace_points = {
    key: generic_curve(K(section[0]) / K(section[2] ** 2), K(section[1]) / K(section[2] ** 3))
    for key, section in trace_sections.items()
}
Q_section = tuple(
    R([
        F(QQ(value).numerator()) / F(QQ(value).denominator())
        for value in pointed["section"][key]
    ])
    for key in (
        "X_coefficients_low_to_high",
        "Y_coefficients_low_to_high",
        "Z_coefficients_low_to_high",
    )
)
Q_generic = generic_curve(
    K(Q_section[0]) / K(Q_section[2] ** 2),
    K(Q_section[1]) / K(Q_section[2] ** 3),
)
M_section = tuple(
    R([
        F(QQ(value).numerator()) / F(QQ(value).denominator())
        for value in bridge_m["section"][key]
    ])
    for key in (
        "X_coefficients_low_to_high",
        "Y_coefficients_low_to_high",
        "Z_coefficients_low_to_high",
    )
)
M_generic = generic_curve(
    K(M_section[0]) / K(M_section[2] ** 2),
    K(M_section[1]) / K(M_section[2] ** 3),
)

expected_trace_profiles = {"S5": (2, 5), "S7": (0, 3), "S17": (0, 3)}
selected_trace_branches = {}
for name, expected in expected_trace_profiles.items():
    matches = [
        branch
        for branch in (1, -1)
        if (
            section_pole_order(trace_sections[(name, branch)]),
            component_depth(trace_sections[(name, branch)]),
        )
        == expected
    ]
    if len(matches) != 1:
        raise ArithmeticError(f"{name} has {len(matches)} branches with profile {expected}")
    selected_trace_branches[name] = matches[0]

s7_point = trace_points[("S7", selected_trace_branches["S7"])]
zero_candidates = []
zero_candidate_audit = []
for trace_signs, q_coefficient in itertools.product(
    itertools.product((1, -1), repeat=3), (-1, 1)
):
    signs = dict(zip(("S5", "S7", "S17"), trace_signs))
    point = (
        signs["S5"] * trace_points[("S5", selected_trace_branches["S5"])]
        - 2 * signs["S7"] * s7_point
        - signs["S17"] * trace_points[("S17", selected_trace_branches["S17"])]
        + q_coefficient * Q_generic
    )
    section = normalized_point(point)
    zero_profile = None if section is None else [
        section_pole_order(section), component_depth(section)
    ]
    difference = normalized_point(point - M_generic) if section is not None else None
    difference_profile = None if difference is None else [
        section_pole_order(difference), component_depth(difference)
    ]
    zero_candidate_audit.append(
        {
            "Q_coefficient_on_stored_positive_section": q_coefficient,
            "trace_signs": signs,
            "zero_profile": zero_profile,
            "Opinned_minus_M_profile": difference_profile,
        }
    )
    if zero_profile == [3, 2] and difference_profile == [5, 3]:
        zero_candidates.append((signs, q_coefficient, point, section, difference))
if len(zero_candidates) != 1:
    raise ArithmeticError(
        f"expected one M-marked pinned zero, got {len(zero_candidates)}; "
        f"audit={zero_candidate_audit}"
    )
selected_trace_signs, q_coefficient, zero_point, zero_section, zero_minus_M = zero_candidates[0]
stored_Q_is_Qminus = q_coefficient == -1
orientation = {
    "S5_branch": selected_trace_branches["S5"],
    "S7_branch": selected_trace_branches["S7"],
    "S17_branch": selected_trace_branches["S17"],
    "Q_coefficient_on_stored_positive_section": q_coefficient,
    "trace_signs": selected_trace_signs,
}
zero_record = {"section": zero_section, "orientations": [orientation]}
selected_signed_trace_sections = {
    name: normalized_point(
        selected_trace_signs[name]
        * trace_points[(name, selected_trace_branches[name])]
    )
    for name in ("S5", "S7", "S17")
}
zero_audit = {
    "individual_trace_profiles": {
        f"{name}:{branch}": [section_pole_order(section), component_depth(section)]
        for (name, branch), section in trace_sections.items()
    },
    "selected_trace_branches": selected_trace_branches,
    "Q_sign_candidate_audit": zero_candidate_audit,
    "stored_Q_is_lattice_Qminus": stored_Q_is_Qminus,
}

zero_inputs = (TRACES, A11, ROUTE, POINTED, BRIDGE_M)
zero_payload = {
    "schema": "elkies-k3.h3-q24-a11-pinned-zero-section-modp.v1",
    "status": "PASS_Q24_A11_PINNED_ZERO_SECTION_RECONSTRUCTION_MODP",
    "prime": int(p),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in zero_inputs],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in zero_inputs
        },
    },
    "lattice_word": route["pinned_zero_low_trace_word"],
    "selected_orientation": orientation,
    "marking_audit": zero_audit,
    "selected_trace_sections": {
        name: section_row(selected_signed_trace_sections[name])
        for name in ("S5", "S7", "S17")
    },
    "section": section_row(zero_section),
    "method": "individual univariate Newton/Pade traces, exact F_p(T) group law, and exact-M sign marking",
    "large_Groebner_required": False,
    "proof_boundary": (
        "Exact over the stated good-reduction field. The exact marked M section "
        "selects between the two Q-sign candidates. Characteristic-zero lifting and literal QQ verification "
        "of O_pinned remain open."
    ),
}
args.zero_output.parent.mkdir(parents=True, exist_ok=True)
args.zero_output.write_text(json.dumps(zero_payload, indent=2, sort_keys=True) + "\n")
print(
    "A11PINNEDZEROMODP|prime={}|degrees={}|PO={}|depth={}|Qminus={}|status={}".format(
        p,
        ",".join(map(str, zero_payload["section"]["degrees_X_Y_Z"])),
        zero_payload["section"]["P_dot_equation_zero"],
        zero_payload["section"]["I12_component_depth_up_to_negation"],
        int(stored_Q_is_Qminus),
        zero_payload["status"],
    ),
    flush=True,
)
print(f"ZERO_OUTPUT|{args.zero_output.resolve()}", flush=True)
if args.stop_after_zero:
    raise SystemExit(0)

zero_generic = zero_point
zero_plus_M = normalized_point(zero_generic + M_generic)
zero_marking_audit = {
    "zero_minus_M_profile": None
    if zero_minus_M is None
    else [section_pole_order(zero_minus_M), component_depth(zero_minus_M)],
    "zero_plus_M_profile": None
    if zero_plus_M is None
    else [section_pole_order(zero_plus_M), component_depth(zero_plus_M)],
}
if (
    zero_minus_M is None
    or section_pole_order(zero_minus_M) != 5
    or component_depth(zero_minus_M) != 3
):
    rejection_inputs = (args.zero_output.resolve(), BRIDGE_M, ROUTE, A11)
    rejection_payload = {
        "schema": "elkies-k3.h3-q24-a11-bridge-m-marking-rejection-modp.v1",
        "status": "REJECT_CURRENT_Q24_A11_BRIDGE_M_MARKING_AT_PINNED_GOOD_REDUCTION",
        "prime": int(p),
        "inputs": {
            "paths": [str(path.relative_to(ROOT)) for path in rejection_inputs],
            "sha256": {
                str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in rejection_inputs
            },
        },
        "claimed_pinned_MW": bridge_m["section"]["pinned_lattice_MW_Abel_Jacobi"],
        "expected_Opinned_minus_M_profile": [5, 3],
        "observed_profiles": zero_marking_audit,
        "stored_Q_is_lattice_Qminus": stored_Q_is_Qminus,
        "exact_modp_group_law": True,
        "large_Groebner_required": False,
        "conclusion": (
            "The rational polynomials remain an exact QQ A11 section, but the "
            "current trace-orientation argument does not mark it as M. Rebuild "
            "the bridge trace with its branches pinned against O_pinned before "
            "using this section in the q8 target construction."
        ),
        "proof_boundary": (
            "Exact section arithmetic at the repository's pinned good prime. "
            "This rejects the claimed pinned-good-reduction MW marking, not the "
            "literal characteristic-zero Weierstrass identity of the section."
        ),
    }
    args.bridge_audit_output.parent.mkdir(parents=True, exist_ok=True)
    args.bridge_audit_output.write_text(
        json.dumps(rejection_payload, indent=2, sort_keys=True) + "\n"
    )
    print(
        "A11BRIDGEMMARKING|expected=5,3|observed={}|status={}".format(
            ",".join(map(str, zero_marking_audit["zero_minus_M_profile"])),
            rejection_payload["status"],
        ),
        flush=True,
    )
    print(f"BRIDGE_AUDIT_OUTPUT|{args.bridge_audit_output.resolve()}", flush=True)
    raise SystemExit(0)

# The trace word S6-2*S2-2*S8 has old-degree weight one.  It therefore
# reconstructs R-a on the C10-zero equation, where R=P12-M; this point has
# pole order five and I12 depth one.  Subtracting O_pinned(eq)=-a recovers the
# genuinely small residual R, with pole order zero and I12 depth three.
residual_samples = {
    orientation: []
    for orientation in itertools.product(
        itertools.product((1, -1), (1, -1)), repeat=3
    )
}
for sample in traces["samples"]:
    tau, curve, points = sample_points(sample)
    for orientation, values in residual_samples.items():
        (branch2, sign2), (branch6, sign6), (branch8, sign8) = orientation
        point = (
            sign6 * points[("S6", branch6)]
            - 2 * sign2 * points[("S2", branch2)]
            - 2 * sign8 * points[("S8", branch8)]
        )
        if point.is_zero():
            continue
        x, y = point.xy()
        values.append((tau, F(x), F(y)))

residual_solutions = {}
residual_audit = []
for orientation, samples in residual_samples.items():
    answers = reconstruct_section(samples, 5)
    accepted = [section for section in answers if component_depth(section) == 1]
    residual_audit.append(
        {
            "orientation": {
                "S2_branch": orientation[0][0],
                "S2_trace_sign": orientation[0][1],
                "S6_branch": orientation[1][0],
                "S6_trace_sign": orientation[1][1],
                "S8_branch": orientation[2][0],
                "S8_trace_sign": orientation[2][1],
            },
            "usable_samples": len(samples),
            "P_dot_O_5_reconstruction_count": len(answers),
            "component_depth_1_count": len(accepted),
        }
    )
    for section in accepted:
        key = tuple(tuple(poly.list()) for poly in section)
        residual_solutions.setdefault(key, {"section": section, "orientations": []})[
            "orientations"
        ].append(orientation)
marked_residuals = []
for record in residual_solutions.values():
    trace_section = record["section"]
    trace_generic = generic_curve(
        K(trace_section[0]) / K(trace_section[2] ** 2),
        K(trace_section[1]) / K(trace_section[2] ** 3),
    )
    section = normalized_point(trace_generic - zero_generic)
    if section is None or [section_pole_order(section), component_depth(section)] != [0, 3]:
        continue
    residual_point = generic_curve(
        K(section[0]) / K(section[2] ** 2),
        K(section[1]) / K(section[2] ** 3),
    )
    target = normalized_point(M_generic + residual_point)
    if target is None or [section_pole_order(target), component_depth(target)] != [11, 4]:
        continue
    target_point = generic_curve(
        K(target[0]) / K(target[2] ** 2), K(target[1]) / K(target[2] ** 3)
    )
    difference = normalized_point(target_point - zero_generic)
    if difference is None or [section_pole_order(difference), component_depth(difference)] != [6, 6]:
        continue
    marked_residuals.append((record, trace_section, section, target, difference))
if len(marked_residuals) != 1:
    raise ArithmeticError(
        f"expected one fully marked P12-M residual, got {len(marked_residuals)} "
        f"from {len(residual_solutions)} trace sections"
    )
(
    residual_record,
    residual_trace_section,
    residual_section,
    target_section,
    target_minus_zero,
) = marked_residuals[0]

payload = {
    "schema": "elkies-k3.h3-q24-a11-q8-horizontal-points-modp.v1",
    "status": "PASS_Q24_A11_Q8_HORIZONTAL_POINTS_RECONSTRUCTION_MODP",
    "prime": int(p),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "equation_zero_translation": route["equation_zero"],
    "q8_divisor": route["q8_divisor"],
    "pinned_zero": {
        "section": section_row(zero_section),
        "selected_orientations": zero_record["orientations"],
        "difference_Opinned_minus_M": section_row(zero_minus_M),
        "orientation_audit": zero_audit,
        "selected_trace_sections": zero_payload["selected_trace_sections"],
    },
    "residual_P12_minus_M": {
        "section": section_row(residual_section),
        "expected_pinned_A11_MW": [-1, 0, -1, 0, 0, 0],
        "degree_weight_one_trace_before_zero_subtraction": section_row(
            residual_trace_section
        ),
        "selected_orientations": [
            {
                "S2_branch": row[0][0],
                "S2_trace_sign": row[0][1],
                "S6_branch": row[1][0],
                "S6_trace_sign": row[1][1],
                "S8_branch": row[2][0],
                "S8_trace_sign": row[2][1],
            }
            for row in residual_record["orientations"]
        ],
        "orientation_audit": residual_audit,
    },
    "q8_target": {
        "section": section_row(target_section),
        "construction": "M(eq)+(P12-M)",
    },
    "marked_difference_P12_minus_Opinned": {
        "section": section_row(target_minus_zero),
        "expected_pinned_A11_MW": route["q8_divisor"]["pinned_target_MW"],
        "P_dot_equation_zero": 6,
        "I12_component_depth_up_to_negation": 6,
    },
    "method": "fibrewise exact trace words followed by univariate Newton/Pade reconstruction and exact mod-p group-law marking",
    "large_Groebner_required": False,
    "proof_boundary": (
        "Exact over the pinned good-reduction field, with exact characteristic-zero "
        "lattice words and zero translation. The displayed horizontal sections "
        "still require characteristic-zero coefficient lifting and literal QQ "
        "verification before they can define the exact q8 H0 plane."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11Q8HORIZONTALMODP|prime={}|Opinned_degrees={}|P12_degrees={}|"
    "residual_degrees={}|difference_degrees={}|status={}".format(
        p,
        ",".join(map(str, payload["pinned_zero"]["section"]["degrees_X_Y_Z"])),
        ",".join(map(str, payload["q8_target"]["section"]["degrees_X_Y_Z"])),
        ",".join(map(str, payload["residual_P12_minus_M"]["section"]["degrees_X_Y_Z"])),
        ",".join(
            map(str, payload["marked_difference_P12_minus_Opinned"]["section"]["degrees_X_Y_Z"])
        ),
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
